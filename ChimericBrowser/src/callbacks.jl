test_callback!(app::Dash.DashApp) = 
callback!(app, Output("info-output", "children"), Input("reads-slider", "value")) do slider_value
    "$slider_value"
end