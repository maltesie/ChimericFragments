test_callback!(app::Dash.DashApp) = 
callback!(app, Input("reads-slider", "value"), Output("info-output", "children")) do slider_value
    "$slider_value"
end